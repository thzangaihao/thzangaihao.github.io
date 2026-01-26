/* ============================================================
   文件名: global.js
   描述: 全站通用脚本 (完美修复版)
   ============================================================ */

/* ------------------------------------------------------------
   1. 全局配置 (Configuration)
   ------------------------------------------------------------ */
window.MathJax = {
    tex: { inlineMath: [['$', '$'], ['\\(', '\\)']], displayMath: [['$$', '$$'], ['\\[', '\\]']] },
    options: { ignoreHtmlClass: 'tex2jax_ignore', processHtmlClass: 'tex2jax_process' }
};

/* ------------------------------------------------------------
   2. 全局工具函数 (Utilities)
   ------------------------------------------------------------ */
function getRelativePath(dbPath) { return '../../' + dbPath; }

let copyFeedbackTimer = null;
let currentActiveElement = null;
let originalContentCache = "";

function resetCopyFeedback() {
    if (currentActiveElement) {
        currentActiveElement.innerHTML = originalContentCache;
        currentActiveElement = null;
        originalContentCache = "";
    }
    if (copyFeedbackTimer) { clearTimeout(copyFeedbackTimer); copyFeedbackTimer = null; }
}

function copyToClipboard(text, element, event) {
    if (event) event.preventDefault();
    if (currentActiveElement && currentActiveElement !== element) resetCopyFeedback();
    if (copyFeedbackTimer) clearTimeout(copyFeedbackTimer);

    navigator.clipboard.writeText(text).then(() => {
        if (currentActiveElement !== element) {
            currentActiveElement = element;
            originalContentCache = element.innerHTML;
        }
        element.innerHTML = '<i class="fas fa-check" style="color: #27ae60;"></i> 已复制';
        copyFeedbackTimer = setTimeout(() => resetCopyFeedback(), 1500);
    }).catch(err => console.error('复制失败:', err));
}

function copyCode(text, btnElement) {
    navigator.clipboard.writeText(text).then(() => {
        const originalHtml = btnElement.innerHTML;
        btnElement.innerHTML = '<i class="fas fa-check"></i>';
        btnElement.classList.add('copied');
        setTimeout(() => {
            btnElement.innerHTML = originalHtml;
            btnElement.classList.remove('copied');
        }, 1500);
    }).catch(err => { console.error('复制失败:', err); btnElement.innerText = 'Error'; });
}

/* ------------------------------------------------------------
   3. 核心 UI 逻辑
   ------------------------------------------------------------ */
function loadFooter() {
    const path = window.footerPath || '/footer.html';
    return fetch(path)
        .then(res => { if (!res.ok) throw new Error(`Footer not found`); return res.text(); })
        .then(data => { const el = document.getElementById('footer-placeholder'); if (el) el.innerHTML = data; })
        .catch(err => console.error('Footer Error:', err));
}

/* ------------------------------------------------------------
   3.2 表格响应式增强 (自动提取标题到外部)
   ------------------------------------------------------------ */
function makeTablesResponsive() {
    document.querySelectorAll('.content table').forEach(table => {
        // 防止重复处理
        if (table.parentElement.classList.contains('table-wrapper')) return;

        // 1. 创建圆角滚动容器
        const wrapper = document.createElement('div');
        wrapper.className = 'table-wrapper';
        
        // 2. [核心升级] 检查并提取表格标题
        const caption = table.querySelector('caption');
        let externalCaption = null;
        
        if (caption) {
            // 创建一个外部的 div 来装标题
            externalCaption = document.createElement('div');
            externalCaption.className = 'table-bottom-caption'; // 新的样式类名
            externalCaption.innerHTML = caption.innerHTML;      // 复制内容
            
            // 从表格里移除原标题 (防止重复显示)
            table.removeChild(caption); 
        }

        // 3. 组装 DOM
        // 先把 wrapper 插到 table 前面
        table.parentNode.insertBefore(wrapper, table);
        // 把 table 移入 wrapper
        wrapper.appendChild(table);
        
        // 4. [核心升级] 把提取出来的标题插到 wrapper 后面
        if (externalCaption) {
            wrapper.parentNode.insertBefore(externalCaption, wrapper.nextSibling);
        }
    });
}

function autoFillArticleInfo() {
    const container = document.querySelector('.article-container');
    if (!container) return;
    const getMeta = (n) => document.querySelector(`meta[name="${n}"]`)?.content;
    const [collId, collTitle, pPath, pTitle, date] = [getMeta('collection'), getMeta('collection-title'), getMeta('parent-path'), getMeta('parent-title'), getMeta('date')];
    if (!collId || !pPath || !date) return;

    const h1 = container.querySelector('h1');
    if (h1) {
        h1.insertAdjacentHTML('beforebegin', `<a href="${pPath}#${collId}" class="back-link"><i class="fas fa-arrow-left"></i> 返回 ${pTitle} / ${collTitle}</a>`);
        h1.insertAdjacentHTML('afterend', `<div class="meta"><span><i class="far fa-calendar"></i> ${date}</span><span style="margin-left:20px"><i class="fas fa-folder"></i> ${collTitle}</span></div>`);
    }
}

/* ------------------------------------------------------------
   4. 导航逻辑
   ------------------------------------------------------------ */
function showCollection(collectionId, btnElement) {
    document.querySelectorAll('.article-list').forEach(el => el.classList.remove('active'));
    const target = document.getElementById(collectionId + '-section');
    if (target) target.classList.add('active');
    document.querySelectorAll('.collection-btn').forEach(el => el.classList.remove('active'));
    if (btnElement) btnElement.classList.add('active');
}

function handleHashNavigation() {
    const hash = window.location.hash.substring(1);
    if (!hash) return;
    const targetBtn = document.querySelector(`button[onclick*="'${hash}'"]`);
    if (targetBtn) showCollection(hash, targetBtn);
}

/* ------------------------------------------------------------
   5. 内容渲染 (卡片与最新文章)
   ------------------------------------------------------------ */

/* [修复] 找回丢失的 renderArticleCards 函数 */
function renderArticleCards() {
    if (!window.ARTICLE_DATABASE || !Array.isArray(window.ARTICLE_DATABASE)) return;
    
    window.ARTICLE_DATABASE.forEach(article => {
        // 寻找带 -section 后缀的容器
        const container = document.getElementById(article.collection + '-section');
        if (!container) return;

        // 生成简单卡片
        const html = `
            <a href="${getRelativePath(article.path)}" class="article-item">
                <div class="article-info">
                    <h4>${article.title}</h4>
                    <p>${article.summary}</p>
                </div>
                <span class="article-date">${article.date}</span>
            </a>`;
        container.insertAdjacentHTML('beforeend', html);
    });
}

/* 渲染首页最新文章 (灰色标签版) */
function renderLatestArticles() {
    const container = document.getElementById('latest-articles-list');
    if (!container) return;

    if (!window.ARTICLE_DATABASE) { container.innerHTML = '<p style="text-align:center;color:#999">暂无数据</p>'; return; }
    
    // 文章显示个数/数目
    const latest = [...window.ARTICLE_DATABASE].sort((a, b) => new Date(b.date) - new Date(a.date)).slice(0, 10);
    
    const html = latest.map(article => {
        let sectionName = '动态'; 
        let iconClass = 'fas fa-newspaper'; 

        if (article.path.includes('dry_lab')) { sectionName = '干实验'; iconClass = 'fas fa-code'; } 
        else if (article.path.includes('wet_lab')) { sectionName = '湿实验'; iconClass = 'fas fa-flask'; } 
        else if (article.path.includes('resources')) { sectionName = '资源站'; iconClass = 'fas fa-book'; } 
        else if (article.path.includes('about')) { sectionName = '关于我'; iconClass = 'fas fa-user'; }

        return `
            <a href="${article.path}" class="article-item">
                <div class="article-info">
                    <span class="article-tag"><i class="${iconClass}"></i> ${sectionName}</span>
                    <h4>${article.title}</h4>
                    <p>${article.summary}</p>
                </div>
                <span class="article-date">${article.date}</span>
            </a>`;
    }).join('');
    container.innerHTML = html;
}

/* ------------------------------------------------------------
   6. 代码高亮与增强
   ------------------------------------------------------------ */
function loadPrismHighlighter() {
    const codeBlocks = document.querySelectorAll('pre code');
    if (codeBlocks.length === 0) return;
    const cdn = 'https://cdnjs.cloudflare.com/ajax/libs/prism/1.29.0';
    const link = document.createElement('link'); link.rel='stylesheet'; link.href=`${cdn}/themes/prism-tomorrow.min.css`; document.head.appendChild(link);
    const script = document.createElement('script'); script.src=`${cdn}/prism.min.js`; 
    script.onload = () => {
        const auto = document.createElement('script'); auto.src=`${cdn}/plugins/autoloader/prism-autoloader.min.js`;
        auto.onload = () => window.Prism.highlightAll();
        document.body.appendChild(auto);
    };
    document.body.appendChild(script);
}

function enhanceCodeBlocks() {
    document.querySelectorAll('pre code').forEach(code => {
        const pre = code.parentElement;
        if (pre.parentElement.classList.contains('code-wrapper')) return;
        const lang = (Array.from(code.classList).find(c => c.startsWith('language-')) || 'TEXT').replace('language-', '').toUpperCase();
        const wrap = document.createElement('div'); wrap.className = 'code-wrapper';
        wrap.innerHTML = `<div class="code-header"><span class="lang-label">${lang}</span><button class="copy-code-btn" onclick="copyCode(this.parentElement.nextElementSibling.innerText, this)"><i class="far fa-copy"></i></button></div>`;
        pre.parentNode.insertBefore(wrap, pre); wrap.appendChild(pre);
    });
}

/* ------------------------------------------------------------
   7. 目录生成 (TOC) - [高级版：防抖 + 智能高亮 + 整体淡入]
   ------------------------------------------------------------ */
function generateTOC() {
    const container = document.querySelector('.article-container');
    const content = document.querySelector('.content');
    if (!container || !content) return;

    const headers = content.querySelectorAll('h2, h3');
    // 如果标题太少，不生成目录，直接显示文章
    if (headers.length < 2) {
        container.style.opacity = '1'; 
        return;
    }

    // --- DOM 构建 ---
    const wrapper = document.createElement('div');
    wrapper.className = 'article-layout-wrapper';
    container.parentNode.insertBefore(wrapper, container);
    
    const sidebar = document.createElement('aside');
    sidebar.className = 'article-sidebar';
    sidebar.innerHTML = `<div class="toc-card"><h3><i class="fas fa-list-ul"></i> 目录</h3><ul id="toc-list"></ul></div>`;
    
    wrapper.appendChild(container);
    wrapper.appendChild(sidebar);

    const list = sidebar.querySelector('#toc-list');
    
    // 侧边栏滚动辅助函数
    const scrollSidebar = (link) => {
        const card = document.querySelector('.toc-card');
        if (!card) return;
        const linkTop = link.offsetTop; 
        const cardHeight = card.clientHeight;
        const scrollTop = card.scrollTop;
        if (linkTop < scrollTop + 20 || linkTop > scrollTop + cardHeight - 20) {
            card.scrollTo({ top: linkTop - 60, behavior: 'smooth' });
        }
    };

    // 生成列表
    headers.forEach((h, i) => {
        if (!h.id) h.id = 'sec-' + i;
        const li = document.createElement('li');
        if (h.tagName === 'H3') li.className = 'toc-sub-item';
        
        const link = document.createElement('a');
        link.href = `#${h.id}`;
        link.className = 'toc-link';
        link.dataset.target = h.id;
        link.textContent = h.textContent;

        // 点击事件 (带锁，防止乱跳)
        link.onclick = (e) => {
            e.preventDefault();
            isClicking = true; 

            document.querySelectorAll('.toc-link').forEach(l => l.classList.remove('active'));
            link.classList.add('active');
            
            scrollSidebar(link);
            document.getElementById(h.id).scrollIntoView({ behavior: 'smooth' });

            if (window.clickTimer) clearTimeout(window.clickTimer);
            window.clickTimer = setTimeout(() => { isClicking = false; }, 1000);
        };

        li.appendChild(link);
        list.appendChild(li);
    });

    // 触发整体淡入动画
    void wrapper.offsetWidth;
    wrapper.classList.add('loaded');

    // --- 滚动监听 (智能版) ---
    let isClicking = false;
    const visibleHeaders = new Set();

    const observer = new IntersectionObserver(entries => {
        // 1. 始终更新可见集合
        entries.forEach(entry => {
            if (entry.isIntersecting) visibleHeaders.add(entry.target.id);
            else visibleHeaders.delete(entry.target.id);
        });

        // 2. 如果正在点击，暂停 UI 更新
        if (isClicking) return;

        // 3. 仲裁：优先选择 DOM 顺序最靠前的可见标题
        let activeId = null;
        for (const header of headers) {
            if (visibleHeaders.has(header.id)) {
                activeId = header.id;
                break; 
            }
        }

        if (activeId) {
            document.querySelectorAll('.toc-link').forEach(l => {
                const isActive = l.dataset.target === activeId;
                l.classList.toggle('active', isActive);
                if (isActive) scrollSidebar(l);
            });
        }
    }, { rootMargin: '0px 0px -80% 0px', threshold: 0 });

    headers.forEach(h => observer.observe(h));
}

/* ------------------------------------------------------------
   8. 初始化入口
   ------------------------------------------------------------ */
window.addEventListener('DOMContentLoaded', async () => {
    if ('scrollRestoration' in history) history.scrollRestoration = 'manual';
    
    await loadFooter();

    // 优先渲染内容
    if (document.body.classList.contains('category-page')) {
        renderArticleCards(); // 现在这个函数存在了，二级页面会正常显示
        handleHashNavigation();
    }
    renderLatestArticles();   
    
    // 装饰性功能
    enhanceCodeBlocks();      
    loadPrismHighlighter(); 
    makeTablesResponsive();
    
    // 文章页逻辑 (TOC 和 信息栏)
    autoFillArticleInfo();
    generateTOC();
});